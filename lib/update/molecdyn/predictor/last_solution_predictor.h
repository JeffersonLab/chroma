// -*- C++ -*-
// $Id: last_solution_predictor.h,v 1.5 2005-04-10 21:46:43 edwards Exp $
/*! \file
 * \brief Last solution predictor
 *
 * Predictors for HMC
 */

#ifndef __last_solution_predictor_h__
#define __last_solution_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace LastSolution4DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  //! Last solution predictor
  /*! @ingroup predictor */
  class LastSolution4DChronoPredictor : 
    public AbsChronologicalPredictor4D<LatticeFermion> {
    
  public:

    // Destructor is automagic
    ~LastSolution4DChronoPredictor(void) {}

    LastSolution4DChronoPredictor(void) : last_solution_available(false) {};

    LastSolution4DChronoPredictor(const LastSolution4DChronoPredictor& p) :
      last_solution(p.last_solution), last_solution_available(p.last_solution_available) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(LatticeFermion& psi,
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi) {
      QDPIO::cout << "LastSolution4DChronoPredictor: ";
      if( last_solution_available ) { 
	QDPIO::cout << "Giving you the last solution" << endl;
	psi = last_solution;
      }
      else {
	QDPIO::cout << "No available last guess. Giving you zero" << endl;
	psi = zero;
      }
    }
    
    // No internal state so reset is a nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;


      // Set the dirty bit
      last_solution_available = false;
    }

    // Ignore new vector
    void newVector(const LatticeFermion& psi) {
      QDPIO::cout << "LastSolutionPredictor: registering new solution" << endl;
      last_solution = psi;
      last_solution_available = true;
    }

  private:

    LatticeFermion last_solution;
    bool last_solution_available;
  };

  

  /*! @ingroup predictor */
  namespace LastSolution5DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };
  
  //! Last solution predictor
  /*! @ingroup predictor */
  class LastSolution5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {

  public:
    ~LastSolution5DChronoPredictor(void) {}

    // Creation
    LastSolution5DChronoPredictor(const int N5_) : N5(N5_), last_solution_available(false) {}

    // Copying
    LastSolution5DChronoPredictor(const LastSolution5DChronoPredictor& p) : 
      N5(p.N5), last_solution(p.last_solution), last_solution_available(p.last_solution_available) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperator< multi1d<LatticeFermion> >& A,
		    const multi1d<LatticeFermion>& chi) { 

      QDPIO::cout << "LastSolutionPredictor:";
      psi.resize(N5);
      if( last_solution_available ) { 
	QDPIO::cout << " last solution is available. Giving you it" << endl;
	psi = last_solution;
      }
      else {
	QDPIO::cout << " last solution is not available. Giving you zero" << endl;
	psi = zero;
      }

    }
    

    // No internal state so reset is a Nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
      last_solution_available = false;
    }

    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) {
      QDPIO::cout << "LastSolutionPredictor: registering new solution" << endl;

      if ( psi.size() != N5 ) { 
	QDPIO::cerr << "Vector of incompatible size presented to Chronological Predictor. Vector.size() = " << psi.size() << " predictor.size()=" << N5 << endl;
	QDP_abort(1);
      }

      last_solution = psi;
      last_solution_available = true;

    }
    
  private:
    const int N5;
    multi1d<LatticeFermion> last_solution;
    bool last_solution_available;

  };
  
}; // End Namespace Chroma

#endif 
