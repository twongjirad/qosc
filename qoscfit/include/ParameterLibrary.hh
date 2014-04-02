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
* \class ParameterLibrary
* \ingroup QoscFit
* \brief Repository for the user to register specific implementation of parameters for use in analyses
*
* -------------------------------------------------------------------------------------------*/

/*
 ------------------------------------------------
 The idea of this class makes me nautious.
 But I do need a way for the python interface
 and other users to add in custom parameters.
 well, we work with compiled code, so the only
 way i can think of to make this work is to
 have a container where predefined parameters
 and the user defined parameters can be stored
 and accessed at run time. so dirty.
 ------------------------------------------------
*/

#ifndef __ParameterLibrary__
#define __ParameterLibrary__

#include "ParameterManager.hh"
#include <string>

namespace qosc {

  class ModelParameter;

  class ParameterLibrary : public ParameterManager {

  public:
  
    static ParameterLibrary* GetTheParameterLibrary();
    static ModelParameter* GetParameter( std::string name );
    static void Print();
    void AddParameter( ModelParameter* parterm, std::string groupname="_nogroup_" );

  protected:

    static ParameterLibrary* singleton;
    ParameterLibrary();
    ParameterLibrary( std::string name );
    virtual ~ParameterLibrary();

    // Parameter Groupings
  protected:
    static std::set< std::string >* m_group_set;
    static std::map< std::string, std::string >* m_groups_dict; // [par,grou]
  public:
    bool IsGroupDefined( std::string groupname );
    std::string GetParameterGroup( std::string parname );
    void GetParametersInGroup( std::string groupname, std::vector< std::string >& parlist );
  };

  ParameterLibrary* getTheParameterLibrary();

}

#endif
