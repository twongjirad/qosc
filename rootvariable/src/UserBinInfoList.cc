#include "UserBinInfoList.hh"
#include <iostream>

using namespace qosc;

UserBinInfoList::UserBinInfoList( std::string binlabel ) {
  m_binlabel = binlabel;
}

UserBinInfoList::~UserBinInfoList() { 
  Clear();
}

UserBinInfoList* UserBinInfoList::Copy( std::string binlabel ) {

  if ( binlabel=="" ) binlabel = GetBinLabel();
  UserBinInfoList* copylist = new UserBinInfoList( binlabel );

  UserBinInfoListIter it;
  for ( it=Begin(); it!=End(); it++) {
    std::string infoname = (*it).first;
    UserBinInfo* info = (*it).second;
    UserBinInfo* copy = info->Copy( binlabel+"_"+infoname );
    //copylist->AddInfo( infoname, copy );
    copylist->AddInfo( infoname, info );
  }

  return copylist;
}

void UserBinInfoList::CopyBinInfoClasses( UserBinInfoList* infolist ) {
  UserBinInfoListIter it;
  for ( it=infolist->Begin(); it!=infolist->End(); it++) {
    std::string infoname = (*it).first;
    UserBinInfo* info = (*it).second;
    UserBinInfo* copy = info->Copy( GetBinLabel()+"_"+infoname );
    AddInfo( infoname, copy );
  }
}

void UserBinInfoList::AddInfo( std::string name, UserBinInfo* info ) { 
  m_infolist[name] = info;
}

void UserBinInfoList::Clear() {
  UserBinInfoListIter it;
  for ( it=Begin(); it!=End(); it++) {
    if ( (*it).second ) delete (*it).second;
  }
  m_infolist.clear();
}

int UserBinInfoList::GetEntries() { 
  return m_infolist.size();
}

void UserBinInfoList::FillInfo( double weight ) {
  if ( GetEntries()<=0 ) return;

  UserBinInfoListIter it;
  for ( it=Begin(); it!=End(); it++) {
    UserBinInfo* info = (*it).second;
    info->FillBinInfo( weight );
  }  
}

UserBinInfo* UserBinInfoList::GetInfo( std::string name ) {
  UserBinInfoListIter it = m_infolist.find( name );
  if ( it==End() ) return NULL;
  return (*it).second;
}

void UserBinInfoList::Write() {
  UserBinInfoListIter it;
  for ( it=Begin(); it!=End(); it++) 
    (*it).second->Write();
}

void UserBinInfoList::Print() {
  std::cout << "Contents of UserBinInfoList[" << GetBinLabel() << "], #instances=" << m_infolist.size() << std::endl;
  for ( UserBinInfoListIter it=Begin(); it!=End(); it++ ) {
    std::cout << "  " << (*it).second->GetInstanceName() << std::endl;
  }
}
