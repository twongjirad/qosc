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

#include "ParameterLibrary.hh"

using namespace qosc;

ParameterLibrary* ParameterLibrary::singleton = NULL;
std::set< std::string >* ParameterLibrary::m_group_set = NULL;
std::map< std::string, std::string >* ParameterLibrary::m_groups_dict = NULL;

ParameterLibrary* getTheParameterLibrary() {
  return ParameterLibrary::GetTheParameterLibrary();
}

ParameterLibrary::ParameterLibrary() 
{}

ParameterLibrary::ParameterLibrary( std::string name ) 
  : ParameterManager( name )
{
  m_group_set = new std::set< std::string >;
  m_groups_dict = new std::map< std::string, std::string >;
  m_group_set->insert( "_nogroup_" );
  m_groups_dict->clear();
}

ParameterLibrary::~ParameterLibrary() {
}

ModelParameter* ParameterLibrary::GetParameter( std::string parname ) {
  ParameterLibrary* parlib = getTheParameterLibrary();
  return dynamic_cast<ParameterManager*>(parlib)->GetParameter( parname );    
}

ParameterLibrary* ParameterLibrary::GetTheParameterLibrary() {
  if ( singleton==NULL )
    singleton = new ParameterLibrary( "ParameterLibrary" );
  return singleton;
}

void ParameterLibrary::AddParameter( ModelParameter* parterm, std::string groupname ) {
  RegisterParameter( parterm );
  m_group_set->insert( groupname );
  std::string parname = parterm->GetName();
  ParameterLibrary::m_groups_dict->insert( std::pair< std::string, std::string>( parname, groupname ) );
}

void ParameterLibrary::Print() {
  ParameterLibrary* parlib = getTheParameterLibrary();
  dynamic_cast<ParameterManager*>(parlib)->Print();
}

bool ParameterLibrary::IsGroupDefined( std::string groupname ) {
  if ( m_group_set->find( groupname )!=m_group_set->end() )
    return true;
  else
    return false;
}


std::string ParameterLibrary::GetParameterGroup( std::string parname ) {
  if ( m_groups_dict->find( parname )!=m_groups_dict->end() )
    return m_groups_dict->find(parname)->second;

  std::cout << "Parameter is not a part of a group." << std::endl;
  assert(false);
}

void ParameterLibrary::GetParametersInGroup( std::string groupname, std::vector< std::string >& parlist ) {
  parlist.clear();
  for ( std::map< std::string, std::string >::iterator igroup=m_groups_dict->begin(); igroup!=m_groups_dict->end(); igroup++ ) {
    if ( igroup->second==groupname )
      parlist.push_back( igroup->first );
  }
}
