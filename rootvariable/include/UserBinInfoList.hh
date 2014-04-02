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
 * \class UserBinInfoList
 * \ingroup rootvariable
 * \brief Container for UserBinInfo objects
 *
 * --------------------------------------------------------------------------------------------
 */

#ifndef __UserBinInfoList__
#define __UserBinInfoList__

#include <map>
#include <string>
#include "UserBinInfo.hh"

namespace qosc {
  // Container for UserBinInfo Objects

  typedef std::map< std::string, UserBinInfo* >::iterator UserBinInfoListIter;

  class UserBinInfoList {

  public:

    UserBinInfoList( std::string binlabel="" );
    ~UserBinInfoList();

    void AddInfo( std::string name, UserBinInfo* info );
    UserBinInfo* GetInfo( std::string name );
    UserBinInfoList* Copy( std::string binlabel="" ); // will create a new container and fill it with deep copies of the user bin info instances stored in it
    void CopyBinInfoClasses( UserBinInfoList* infolist ); // will copy the objects from the BinInfoList instance passed to it
    UserBinInfoListIter Begin() { return m_infolist.begin(); };
    UserBinInfoListIter End() { return m_infolist.end(); };
    void Clear();
    int GetEntries();
    void FillInfo( double weight = 1.0 );
    void Write();
    void SetBinLabel( std::string binlabel ) { m_binlabel = binlabel; };
    std::string GetBinLabel() { return m_binlabel; };
    void Print();
  
  protected:

    std::string m_binlabel; /// this is a way to give this bin info list a label to help organize the instances
    std::map< std::string, UserBinInfo* > m_infolist; /// dictionary of user info instances
  

  };
}
#endif

