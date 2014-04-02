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
 * \class UserBinInfo
 * \ingroup rootvariable
 * \brief Abstract base class used to store user info for each bin in a Hist object.
 *
 * This class is responsible for reweighting a bin in the Hist object.
 * In order to do this, the class is suppose to store the information it needs to reweight each bin.
 * The Hist object will then call for each bin the GetBinReweight function and 
 *   multiply its contents by that returned weight.
 *  
 * --------------------------------------------------------------------------------------------
 */

#ifndef __UserBinInfo__
#define __UserBinInfo__

#include <string>
class TFile;

namespace qosc {
  class UserBinInfo {

  public:

    UserBinInfo( std::string tag="" );
    virtual ~UserBinInfo();

    virtual double GetBinReweight() = 0; /// user must define this
    virtual void FillBinInfo( double weight=1.0 ) = 0;  // list of root variable quantities that the info needs to do its job
    virtual UserBinInfo* Copy( std::string tag="" ) = 0; //  create copy of the bin info class!
    virtual void Write() {}; // provides method for storage to TFile
    virtual void LoadBinInfoFromFile( TFile* file ) {}; // loads information from TFile back into Instance of BinInfo object
    virtual void Print(); // dump information
    virtual void Reset() {}; // reset the instance

  protected:
    std::string m_instance_name; ///< unique specifier.
  public:
    std::string GetInstanceName() { return m_instance_name; };
    std::string GetInstanceNameNoIDtag();

  private:
    int m_instanceID;
    static int __nInstances__; // note the underscores!! don't mess with this.

  public:
    int GetInstanceID() { return m_instanceID; };

  };
}
#endif
