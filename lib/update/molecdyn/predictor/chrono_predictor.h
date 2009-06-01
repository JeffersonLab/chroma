// -*- C++ -*-
// $Id: chrono_predictor.h,v 3.2 2009-06-01 20:39:36 bjoo Exp $
/*! \file
 * \brief Chronological predictor for HMC
 *
 * Chronological predictor for HMC
 */

#ifndef chrono_predictor_h
#define chrono_predictor_h

#include "chromabase.h"
#include "linearop.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma 
{
   
  //! Abstract interface for a Chronological Solution predictor
  /*! @ingroup predictor */
  template<typename T>
  class AbsChronologicalPredictor4D {
  public:
    
    // Virtual destructor to help with cleanup
    virtual ~AbsChronologicalPredictor4D(void) {}

    // Set psi to be the next initial guess
    //
    // I have expanded the interface to allow us to 
    // pass the Matrix M, and the RHS chi as well as phi
    // 
    // We are trying to solve the system: 
    //            A  psi = chi
    //
    // for a CG situation A = MdagM
    //
    // and we are trying to get a guess for phi which 
    // minimises the initial residual.
    virtual void operator()(T& psi, 
			    const LinearOperator<T>& A, 
			    const T& chi) = 0;

    // Reset internal state (call this if the gauge field or 
    // pseudofermion fields change)
    virtual void reset(void) = 0;

    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const T& psi) = 0;
  };

  //! Abstract interface for a Chronological Solution predictor
  /*! @ingroup predictor */
  template<typename T>
  class AbsTwoStepChronologicalPredictor4D : public AbsChronologicalPredictor4D<T> {
  public:
    
    // Virtual destructor to help with cleanup
    virtual ~AbsTwoStepChronologicalPredictor4D(void) {}

    // Set psi to be the next initial guess
    //
    // I have expanded the interface to allow us to 
    // pass the Matrix M, and the RHS chi as well as phi
    // 
    // We are trying to solve the system: 
    //            A  psi = chi
    //
    // for a CG situation A = MdagM
    //
    // and we are trying to get a guess for phi which 
    // minimises the initial residual.


    
    virtual void predictX(T& X, 
			  const LinearOperator<T>& A, 
			  const T& chi) = 0;


    virtual void predictY(T& Y, 
			  const LinearOperator<T>& A, 
			  const T& chi) = 0;

    // Reset internal state (call this if the gauge field or 
    // pseudofermion fields change)
    virtual void reset(void) = 0;


    virtual void newXVector(const T& X) = 0;
    virtual void newYVector(const T& Y) = 0;


    // This is a 'predict X' which is always MdagM X = phi
    // These two are backward compatibilities until everything 
    // is changed...
    virtual void operator()(T& psi, 
			    const LinearOperator<T>& A, 
			    const T& chi) 
    {

      predictX(psi, A, chi);
    }

    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const T& psi) 
    {
      newXVector(psi);
    }

    // These two deal with the 
  };
 

  //! Abstract interface for a Chronological Solution predictor in 5D
  /*! @ingroup predictor */
  template<typename T>
  class AbsChronologicalPredictor5D {
  public:

    virtual ~AbsChronologicalPredictor5D(void) {}

    // Set psi to be the next initial guess
    // I have expanded the interface to allow us to 
    // pass the Matrix M, and the RHS chi as well as phi
    // 
    // We are trying to solve the system
    //            A psi = chi
    // 
    // (for now in a CG based situation A = MdagM)
    //
    // and we are trying to get a guess for phi which 
    // minimises the initial residual.
    virtual void operator()(multi1d<T>& psi,
			    const LinearOperatorArray<T>& A, 
			    const multi1d<T>& chi) = 0;
			    

    // Reset internal state (call this if the gauge field or 
    // pseudofermion fields change)
    virtual void reset(void) = 0;
    
    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const multi1d<T>& psi) = 0;
  };

}; // End namespace
#endif
