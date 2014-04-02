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
