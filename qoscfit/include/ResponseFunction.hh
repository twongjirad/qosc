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
 * \class ResponseFunction
 * \ingroup QoscFit
 * \brief Pure virtual class that defines interface to ResponseFunction objects
 *
 * -------------------------------------------------------------------------------------------*/

#ifndef __ResponseFunction__
#define __ResponseFunction__

#include <string>

namespace qosc {

  class ResponseFunction {

  public:

    ResponseFunction();
    virtual ~ResponseFunction();

    virtual double GetResponse( double x[] ) = 0;
    virtual double GetDerivative( double x[] ) = 0;
    virtual void Add( ResponseFunction* addition, double weight ) = 0;
    virtual std::string IsA() = 0;

  };

}

#endif
