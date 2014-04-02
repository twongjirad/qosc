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
 * -------------------------------------------------------------
 * \class Weight
 * \ingroup rootvariabl
 * \brief This abstract base class represents a generic weight.
 * 
 * Doesn't care what inputs are needed.  Derived classes must
 * implement those details.  But all objects which use a weight
 * needs to be able to calculate a value.
 *
 * --------------------------------------------------------------
 */

#ifndef __Weight__
#define __Weight__

class TChain;

namespace qosc {

  class Weight {

  public:

    Weight();
    virtual ~Weight();

    virtual double CalculateWeight() = 0;
    virtual double CalculateWeight( double variable ) = 0;

    virtual Weight* CloneWeight( TChain* source );
  };

}

#endif

