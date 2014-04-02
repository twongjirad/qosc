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
 * -----------------------------------------------------------------------------------------------------------------------
 * \class OscTerm
 * \defgroup T2KOscFit
 * \brief Model Parameter responsible for oscillation fitting
 *
 * Right now just stores value. Not connected to Weighting interface.  Have to think about this.
 *
 * -----------------------------------------------------------------------------------------------------------------------
 */

#ifndef __OscTerm__
#define __OscTerm__

#include "BasicParameter.hh"

namespace qosc {

  class OscTerm : public BasicParameter {

  public:
    OscTerm( std::string parName, double initial, double lowbound, double upbound );
    virtual ~OscTerm();

    virtual std::string IsA() { return "OscTerm"; };

    virtual void SetValue( double value );  
  };

}

#endif
