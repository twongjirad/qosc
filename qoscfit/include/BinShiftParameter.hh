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
/** --------------------------------------------------------------------------------------------
 * \class BinShiftParameter
 * \ingroup QoscFit
 * \brief Parameter for 1D bin migration
 * 
 * This class represent parameters which
 * perform a 1D bin migrations. It is a MINUIT 
 * term by definition.
 * 
 * ------------------------------------------------------------------------------------------ */

#ifndef __BinShiftParameter__
#define __BinShiftParameter__

#include <string>
#include "BasicParameter.hh"

class TH1D;

namespace qosc {

  class BinShiftParameter : public BasicParameter {

  public:  
    BinShiftParameter( std::string name, std::string sample_to_apply, double mean, double sig );
    virtual ~BinShiftParameter();

    virtual std::string IsA() { return "BinShiftParameter"; };
    virtual void TransformHistogram( TH1D* h );
    virtual void AddAffectedSample( std::string name );

  protected:
    std::set< std::string > m_sample_applies;
  public:
    bool DoesSampleApply( std::string sample );

  protected:
    TH1D* output_hist;

  };

}

#endif
