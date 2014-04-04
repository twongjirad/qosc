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
/* --------------------------------------------------------------------------------------------
 * \class ResponseFunctionPoly
 * \ingroup QoscFit
 * \brief Polynomial response function
 *
 * -------------------------------------------------------------------------------------------*/

#ifndef __ResponseFunctionPoly__
#define __ResponseFunctionPoly__

#include "ResponseFunction.hh"

namespace qosc {

  class ResponseFunctionPoly : public ResponseFunction { 

  public:

    ResponseFunctionPoly( double pars[], int nterms );
    virtual ~ResponseFunctionPoly();

    virtual double GetResponse( double x[] );
    virtual double GetDerivative( double x[] );
    virtual void Add( ResponseFunction* addition, double weight );
    virtual std::string IsA() { return "Poly"; };

    int m_nterms;
    double* params;
    double m_weight;

  };

}

#endif
