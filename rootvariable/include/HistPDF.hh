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
 * \class HistPDF
 * \ingroup rootvariable
 * \brief Interface to PDFBinned objects where a probability needs to be generated from the histogram.
 *
 * Adds the GetProbability() function to the PDFBinned object.
 *
 * --------------------------------------------------------------------------------------------
 */

#ifndef __HistPDF__
#define __HistPDF__

#include "Hist.hh"

class TH1D;
class TH2D;

namespace qosc {

  class HistPDF : public Hist {

  public:
  
    HistPDF( std::string name );
    HistPDF( std::string name, TH1D* hist );
    HistPDF( std::string name, TH2D* hist );
    virtual ~HistPDF();

    virtual double GetProbability( double value );

  protected:

  };

}

#endif
