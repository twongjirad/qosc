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

#include "ResponseFunctionSpline.hh"
#include <assert.h>

#include "TSpline.h"

#include "ResponseFunctionPoly.hh" // this is ugly how the different types of splines need to know about one another

using namespace qosc;

int ResponseFunctionSpline::ninstances = 0;

ResponseFunctionSpline::ResponseFunctionSpline( TSpline3* spline ) {
  m_spline = spline;
  fSplineReady = true;
  m_weight = 0.;
  nknots = m_spline->GetNp();
  x = NULL; y=NULL;
  m_id = ResponseFunctionSpline::ninstances;
  ResponseFunctionSpline::ninstances++;
}

ResponseFunctionSpline::ResponseFunctionSpline() {
  // The empty constructor. Builds an empty spline we can use later.
  m_spline = NULL;
  fSplineReady = false;
  m_weight = 0.;
  nknots = 0;
  x = NULL; y=NULL;
  m_id = ResponseFunctionSpline::ninstances;
  ResponseFunctionSpline::ninstances++;
}

ResponseFunctionSpline::~ResponseFunctionSpline() {
  if ( m_spline ) delete m_spline;
  if ( x ) delete x;
  if ( y ) delete y;
}


void ResponseFunctionSpline::ReadySpline() {
  if ( fSplineReady ) return;
  if ( !x || !y ) assert(false);
  char title[50];
  sprintf( title, "response_spline_%d", m_id );
  m_spline = new TSpline3( title, x, y, nknots );
  fSplineReady = true;
  // don't need the arrays any more.
  delete x;
  delete y;
  x = NULL;
  y = NULL;
}

double ResponseFunctionSpline::GetResponse( double x[] ) {
  if ( !fSplineReady ) ReadySpline();
  return m_spline->Eval( x[0] );
}

double ResponseFunctionSpline::GetDerivative( double x[] ) {
  if ( !fSplineReady ) ReadySpline();
  return m_spline->Derivative( x[0] );
}

void ResponseFunctionSpline::CopyStructure( TSpline3* spline ) {
  nknots = spline->GetNp();
  x = new double[nknots];
  y = new double[nknots];
  for (int i=0; i<nknots; i++) {
    double ix, iy;
    spline->GetKnot( i, ix, iy );
    x[i] = ix;
    y[i] = 0.;
  }
  m_weight = 0.;
}

// ----------------------------------------------------
// Combining Splines and Poly with on another. Other option is to define new class as combination
void ResponseFunctionSpline::Add( ResponseFunction* addition, double weight ) {
  if ( addition->IsA()=="Spline" ) {
    ResponseFunctionSpline* sresponse = dynamic_cast< ResponseFunctionSpline* >( addition );
    if ( !sresponse ) assert(false);
    if ( x==NULL ) CopyStructure( sresponse->m_spline );
    if ( sresponse->m_spline->GetNp()!=nknots ) assert(false);
    double totalw = weight + m_weight;
    if ( totalw>0. ) {
      for ( int i=0; i<nknots; i++ ) {
	double ix, iy;
	sresponse->m_spline->GetKnot( i, ix, iy );
	y[i] = (m_weight*y[i] + weight*iy)/totalw;
      }
      m_weight = totalw;
    }
  }
  else if ( addition->IsA()=="Poly" ) {
    ResponseFunctionPoly* sresponse = dynamic_cast< ResponseFunctionPoly* >( addition );
    double totalw = weight + m_weight;
    if ( totalw>0. ) {
      for ( int i=0; i<nknots; i++ ) {
        double ix, iy;
        //m_spline->GetKnot( i, ix, iy );
	ix = x[i];
	double inputx[1] = { ix };
	iy = sresponse->GetResponse( inputx );
        y[i] = (m_weight*y[i] + weight*iy)/totalw;
      }
      m_weight = totalw;
    }
  }
  else
    assert(false);
}
