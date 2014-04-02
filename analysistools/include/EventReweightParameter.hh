//-*- mode:c++; c-basic-offset:2;   -*-
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
/** ----------------------------------------------------------------------------------
* \class EventReweightParameter
* \defgroup analysistools
*
* Concrete implementation of ModelParameter class.
* This class represents those types of model parameters which modify a MC expectation
*  by reweighting the MC events event-by-event
* The class assumes that there can be more than one 'term' which modifies the reweighting
*
* This class modifies the expectation via reweights to each MC event
* Consequently, the parameter will need it's own weight generator.
* Because the same parameter can, in principle, have a different effect for each type
*  of sample bins, we store a generator for each sample.
* The generators are stored in std::map< std::string, SysTermI* > m_systermI_dict
*
* The SysTermI class is used to interface between this class and some other
*  event reweighting package. Please refer to SysTermI.hh to see what
*  functions are expected to be implemented
*
* ----------------------------------------------------------------------------------*/

#ifndef __EventReweightParameter__
#define __EventReweightParameter__

#include "SysTermI.hh"
#include "BasicParameter.hh"

namespace qosc {

  class EventReweightParameter : public BasicParameter {

  public:
    EventReweightParameter( std::string systermname, double mean, double sig );
    virtual ~EventReweightParameter();

    virtual std::string IsA() { return "EventReweightParameter"; };

    // Sys Term Sample Management
  protected:
    std::set< std::string > m_sample_list; //< Analysis samples sys term applies to
  public:
    void AddSampleTermAppliesTo( std::string sample ) {
      m_sample_list.insert( sample );
    };
    bool DoesTermApplyToSample( std::string sample ) { 
      std::set<std::string>::iterator it = m_sample_list.find( sample );
      if ( it!=m_sample_list.end() ) return true;
      else return false;
    };

    // Event weight generators. One for each sample.
  protected:
    std::map< std::string, SysTermI*> m_systermI_dict; //< <samplename, SysTerm Calculator>: dict to classes that calculate weight of event due to sys term shift: 
    typedef std::map< std::string, SysTermI*>::iterator SysTermIDictIter;
    SysTermIDictIter SysTermIDictBegin() { return m_systermI_dict.begin(); };
    SysTermIDictIter SysTermIDictEnd() { return m_systermI_dict.end(); };

  public:
    virtual void SetEventWeightGenForSample( std::string sample, SysTermI* systerm ) { 
      AddSampleTermAppliesTo( sample );    
      m_systermI_dict[sample] = systerm; 
    };
    SysTermI* GetEventWeightGen( std::string samplename ) {
      if ( m_systermI_dict.find( samplename )!=m_systermI_dict.end() ) return m_systermI_dict[samplename];
      else return NULL;
    };
    bool HasEventWeightGenBeenDefinedForSample( std::string sample ) {
      if ( GetEventWeightGen( sample ) ) return true;
      return false;
    };

    // Overload SetPullValue to pass pull value to event generators
  public:
    virtual void SetPullValue( double pull );

  };
}

#endif
