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

#include "UserBinInfo.hh"
#include <iostream>
#include <sstream>

using namespace qosc;

int UserBinInfo::__nInstances__ = 0;

UserBinInfo::UserBinInfo( std::string tag ) { 
  m_instanceID = __nInstances__;
  __nInstances__++;
  
  std::stringstream ss;
  ss.str("");
  ss << tag << "__xx" << GetInstanceID() << "xx";
  m_instance_name = ss.str(); 

}

UserBinInfo::~UserBinInfo() {
}

void UserBinInfo::Print() {
  std::cout << "UserBinInfo[" << m_instance_name << "] " << this << std::endl;
}

std::string UserBinInfo::GetInstanceNameNoIDtag() {
  size_t pos = m_instance_name.rfind("__xx");
  if ( pos==std::string::npos ) {
    std::cout << "UserBinInfo::GetInstanceNameNoIDtag() -- missing instance tag!!" << std::endl;
    return m_instance_name;
  }
  return m_instance_name.substr(0,pos);
}
