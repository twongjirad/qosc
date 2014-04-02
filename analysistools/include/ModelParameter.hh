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
* \class ModelParameter
* \defgroup AnalysisTools
*
* We expand the idea of this object to be any parameter that can modify the MC
*   expectation.
* 
* Now, we expect that a user will have to define a whole set of different types of 
*   parameters, each with it's own information, and behavior for modifying a MC expectation.
* Because of this, this is class is left as an abstract virtual base class.
* The required function the user must define is:
*   std::string ModelParameter::IsA()
* This is suppose to return the name for the type of concrete systerm object.  
* Also, when the MC expectations (represented by the Sample class type and stored in the SampleManager)
*   are called to be modified, we provide no set rule to do so, but merely pass the container
*   for the sys terms. The user is left to incoporate the information in the systerm class as she 
*   sees fit.
* So what does this class do?
* Well, the AnalysisBase class needs to be able to track and organize systematic terms. It also needs to be 
*   able to sync the fitter and its parameter values to that of this list. So the bulk of the code found here
*   and in this class's container is just for book keeping and methods for the AnalysisBase class.
* So the member information stored here:
*   Type name: get via IsA();
*   Token name: get via GetName();
*   Parameter value: GetValue();
* 
*
// ---------------------------------------------------------------------------------*/


#ifndef __ModelParameter__
#define __ModelParameter__

#include <string>
#include <map>
#include <set>
#include <cmath>
#include <iostream>

namespace qosc {

  typedef enum FitterParType { kNewtonTerm, kMinuitTerm }; 

  class ModelParameter {


  public:

    ModelParameter();
    ModelParameter( std::string systermName,  double sigma, FitterParType t);
    virtual ~ModelParameter();

    virtual std::string IsA() = 0; ///< This is to help user identify among different concrete child types

    // Parameter Token Name
  protected:
    std::string m_parameterName;
  public:
    std::string GetName() { return m_parameterName; };

    // Parameter Fitter Type
  protected:
    FitterParType fTermType;
    void SetFitterParType( FitterParType t ) { fTermType = t; };
  public:
    FitterParType GetFitterParType() { return fTermType; };

    // Parameter Value
  protected:
    double m_parameter_value; ///< pull value given to parameter.
  public:
    virtual double GetValue() { return m_parameter_value; };
    virtual void SetValue( double value ); ///< virtual so that side effects to setting the pull value can be implemented. be careful, though.
    double* GetParameterAddress() { return &m_parameter_value; }; ///< makes me queasy. don't be evil.

    // Parameter Sigma
  protected:
    double m_sigma; ///< size of sigma in units of variable
  public:
    virtual void SetSigma( double sig ) { m_sigma = sig; };
    virtual double GetSigma() { return m_sigma; };

    // Parameter Central Value (or default value)
  protected:
    double m_central_value;
  public:
    virtual void SetCentralValue( double mean ) { m_central_value = mean; };
    virtual double GetCentralValue() { return m_central_value; };

    // Parameter Bounds and Scaling Factor
  protected:
    bool fBounded;
    double lowBound;
    double highBound;
    double m_length_scale;
  public:
    void SetBounds( double low, double high, double length_scale=-1.0 );
    bool IsBounded() { return fBounded; };
    void Unbound() { fBounded = false; };
    double GetHighBound() { return highBound; };
    double GetLowBound() { return lowBound; };
    double GetLengthScale() { return m_length_scale; };
    void SetLengthScale( double length_scale ) { m_length_scale = length_scale; };

    // Parameter Active Status
  protected:
    bool m_active_flag;
  public:
    void SetActiveFlag( bool active ) { m_active_flag = active; };
    bool IsActive() { return m_active_flag; };

    // Parameter Float/Fix status
  protected:
    bool fVaryPar;
  public:
    bool IsVaried() { return fVaryPar; };
    void VaryParameter() { fVaryPar = true; };
    void FixParameter() { fVaryPar = false; };
  
  
    // Parameter ID number. For book keeping by the fitter class.
  protected:
    int m_parameter_id; //< for book keeping. Unindexed (sentinel) value is -1.
  public:
    void SetID( int id ) { m_parameter_id = id; };
    int GetID() { return m_parameter_id; };

    // Parameter Rescaling Factor: this is to help make parameters similar in scale
  protected:
    double m_ScalingFactor;
  public:
    double GetScalingFactor() { return m_ScalingFactor; };
    void SetScalingFactor( double scale ) { 
      // rescale variable
      double rescale = m_ScalingFactor/scale;
      SetBounds( GetLowBound()*rescale, GetHighBound()*rescale );
      SetValue( GetValue()*rescale );
      SetSigma( GetSigma()*rescale );
      SetCentralValue( GetCentralValue()*rescale );
      // set scale
      m_ScalingFactor = scale;
    };

    double GetUnscaledValue() {
      return GetValue()*m_ScalingFactor;
    }
    void SetUnscaledValue( double value ) {
      value /= m_ScalingFactor;
      SetValue( value );
    };


    // Optional default value
  protected:
    double m_default_val;
    bool has_default;
  public:
    void SetDefault(double v);
    double GetDefault();
    void SetToDefault() {
      SetValue( GetDefault() );
    };

    // ------------------------------------------------------------------------------------------------
    // behavior flags
  protected:
    bool fDoesTermPull;
  public:
    void TermDoesNotPull() { fDoesTermPull = false; };
    bool DoesTermPull() { return fDoesTermPull; };

  protected:
    bool fDoesTermThrow;
  public:
    void TermDoesNotThrow() { fDoesTermPull = false; };
    void TermDoesThrow() { fDoesTermPull = true; };
    bool DoesTermThrow() { return fDoesTermPull; };
  


  };

}

#endif
