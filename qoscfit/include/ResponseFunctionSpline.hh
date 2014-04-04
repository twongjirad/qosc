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
/**
 * --------------------------------------------------------------------------------------------
 * \class ResponseFunctionSpline
 * \ingroup QoscFit
 * \brief Concrete implementation of a parameter using spline function
 *
 * -------------------------------------------------------------------------------------------*/

#ifndef __ResponseFunctionSpline__
#define __ResponseFunctionSpline__

#include "ResponseFunction.hh"

class TSpline3;

namespace qosc {

  class ResponseFunctionSpline : public ResponseFunction {

  public:
    ResponseFunctionSpline( TSpline3* spline );
    ResponseFunctionSpline();
    virtual ~ResponseFunctionSpline();

    virtual double GetResponse( double x[] );
    virtual double GetDerivative( double x[] );
    virtual void Add( ResponseFunction* addition, double weight );
    virtual std::string IsA() { return "Spline"; };

    void CopyStructure( TSpline3* spline );
    void ReadySpline();

    TSpline3* m_spline;
    bool fSplineReady;
    // For spline building
    int nknots;
    double m_weight;
    double* x;
    double* y;

  protected:
    int m_id;
    static int ninstances;
  };
}

#endif
