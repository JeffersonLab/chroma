// -*- C++ -*-
// $Id: last_solution_predictor.h,v 3.5 2009-06-01 20:39:37 bjoo Exp $
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
  namespace LastSolution4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Last solution predictor
  /*! @ingroup predictor */
  class LastSolution4DChronoPredictor : 
    public AbsTwoStepChronologicalPredictor4D<LatticeFermion> 
  {
   
  public:

    // Destructor is automagic
    ~LastSolution4DChronoPredictor(void) {}

    LastSolution4DChronoPredictor(void) : last_solutionX_available(false),
					  last_solutionY_available(false),
					  last_solutionX(zero),
					  last_solutionY(zero){}

    LastSolution4DChronoPredictor(const LastSolution4DChronoPredictor& p) :
      last_solutionX_available(p.last_solutionX_available),
      last_solutionY_available(p.last_solutionY_available),
      last_solutionX(p.last_solutionX), 
      last_solutionY(p.last_solutionY) {}

    // Zero out psi -- it is a zero guess after all
    void predictX(LatticeFermion& X,
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi) 
    {
      START_CODE();

      QDPIO::cout << "LastSolution4DChronoPredictor: ";
      if( last_solutionX_available ) { 
	QDPIO::cout << "Giving you the last solution" << endl;
	X = last_solutionX;
      }
      else {
	QDPIO::cout << "No available last guess. Giving you zero" << endl;
	X = zero;
      }
    
      END_CODE();
    }

    // Zero out psi -- it is a zero guess after all
    void predictY(LatticeFermion& Y,
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi) 
    {
      START_CODE();

      QDPIO::cout << "LastSolution4DChronoPredictor: ";
      if( last_solutionY_available ) { 
	QDPIO::cout << "Giving you the last solution" << endl;
	Y = last_solutionY;
      }
      else {
	QDPIO::cout << "No available last guess. Giving you zero" << endl;
	Y = zero;
      }
    
      END_CODE();
    }
    
    // No internal state so reset is a nop
    void reset(void) 
    {

      // Set the dirty bit
      last_solutionX_available = false;
      last_solutionY_available = false;

    }

    // Ignore new vector
    void newXVector(const LatticeFermion& X) 
    {
      START_CODE();

      QDPIO::cout << "LastSolutionPredictor: registering new solution" << endl;
      last_solutionX = X;
      last_solutionX_available = true;
    
      END_CODE();
    }

    // Ignore new vector
    void newYVector(const LatticeFermion& Y) 
    {
      START_CODE();

      QDPIO::cout << "LastSolutionPredictor: registering new solution" << endl;
      last_solutionY = Y;
      last_solutionY_available = true;
    
      END_CODE();
    }

  private:

    bool last_solutionX_available;
    bool last_solutionY_available;
    LatticeFermion last_solutionX;
    LatticeFermion last_solutionY;

  };

  

  /*! @ingroup predictor */
  namespace LastSolution5DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }
  
  //! Last solution predictor
  /*! @ingroup predictor */
  class LastSolution5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> 
  {

  public:
    ~LastSolution5DChronoPredictor(void) {}

    // Creation
    LastSolution5DChronoPredictor(const int N5_) : N5(N5_), last_solution_available(false) {}

    // Copying
    LastSolution5DChronoPredictor(const LastSolution5DChronoPredictor& p) : 
      N5(p.N5), last_solution(p.last_solution), last_solution_available(p.last_solution_available) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& A,
		    const multi1d<LatticeFermion>& chi)
    { 
      START_CODE();

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

      END_CODE();
    }
    

    // No internal state so reset is a Nop
    void reset(void) {
      last_solution_available = false;
    }

    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) 
    {
      START_CODE();

      QDPIO::cout << "LastSolutionPredictor: registering new solution" << endl;

      if ( psi.size() != N5 ) { 
	QDPIO::cerr << "Vector of incompatible size presented to Chronological Predictor. Vector.size() = " << psi.size() << " predictor.size()=" << N5 << endl;
	QDP_abort(1);
      }

      last_solution = psi;
      last_solution_available = true;

      END_CODE();
    }
    
  private:
    const int N5;
    multi1d<LatticeFermion> last_solution;
    bool last_solution_available;

  };
  
} // End Namespace Chroma

#endif 
