// -*- C++ -*-
// $Id: linear_extrap_predictor.h,v 3.5 2009-06-01 20:39:37 bjoo Exp $
/*! \file
 * \brief Linear extrapolation predictor
 *
 * Predictors for HMC
 */

#ifndef __linear_extrap_predictor_h__
#define __linear_extrap_predictor_h__

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/circular_buffer.h"

namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace LinearExtrapolation4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Last solution predictor
  /*! @ingroup predictor */
  class LinearExtrapolation4DChronoPredictor : 
    public AbsTwoStepChronologicalPredictor4D<LatticeFermion> 
  {
  private:
    Handle< CircularBuffer<LatticeFermion> > chrono_bufX;
    Handle< CircularBuffer<LatticeFermion> > chrono_bufY;

  public:
    LinearExtrapolation4DChronoPredictor(void) : chrono_bufX(new CircularBuffer<LatticeFermion>((unsigned int)2)), chrono_bufY(new CircularBuffer<LatticeFermion>((unsigned int)2)) {}

    // Destructor is automagic
    ~LinearExtrapolation4DChronoPredictor(void) {}


    void guess(LatticeFermion& psi,
	       const LinearOperator<LatticeFermion>& A, 
	       const LatticeFermion& chi, 
	       Handle< CircularBuffer <LatticeFermion> > chrono_buf ) 
    {
      START_CODE();

      switch( chrono_buf->size() ) { 
      case 0:
      {
	QDPIO::cout << "LinearExtrapolationPredictor: giving you zero" << endl;
	psi = zero;
      }
      break;

      case 1:
      { 
	QDPIO::cout << "LinearExtrapolationPredictor: giving you last soln" << endl;
	chrono_buf->get(0,psi);
      }
      break;

      case 2:
      {
	QDPIO::cout << "LinearExtrapolationPredictor: giving you linear extrapolation" << endl;

	LatticeFermion y0; 
	chrono_buf->get(0,y0); // Most recent

	LatticeFermion y1; 
	chrono_buf->get(1,y1); // Least recent

	psi = Real(2)*y0 - y1;         // Linear Extrapolation

      }
      break;
      default:
	QDPIO::cerr << "Unknown case reached in LinearExtrapPredictor " << endl;
	QDP_abort(1);
	break;
      }

      END_CODE();
    }

    void predictX(LatticeFermion& X, 
		  const LinearOperator<LatticeFermion>& A, 
		  const LatticeFermion& chi) {
      guess(X,A,chi,chrono_bufX);
    }

    void predictY(LatticeFermion& Y, 
		  const LinearOperator<LatticeFermion>& A, 
		  const LatticeFermion& chi) {
      guess(Y,A,chi,chrono_bufY);
    }

    
    // No internal state so reset is a nop
    void reset(void) {
      chrono_bufX->reset();
      chrono_bufY->reset();
    }

    // Ignore new vector
    void newXVector(const LatticeFermion& X) 
    {
      START_CODE();

      QDPIO::cout << "LinearExtrapolationPredictor: registering new X solution. ";
      chrono_bufX->push(X);
      QDPIO::cout << " number of vectors stored is = " << chrono_bufX->size() << endl;
    
      END_CODE();
    }

    // Ignore new vector
    void newYVector(const LatticeFermion& Y) 
    {
      START_CODE();

      QDPIO::cout << "LinearExtrapolationPredictor: registering new Y solution. ";
      chrono_bufY->push(Y);
      QDPIO::cout << " number of vectors stored is = " << chrono_bufY->size() << endl;
    
      END_CODE();
    }


  };

  

  /*! @ingroup predictor */
  namespace LinearExtrapolation5DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }
  
  //! Last solution predictor
  /*! @ingroup predictor */
  class LinearExtrapolation5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> 
  {
  private: 
    Handle< CircularBufferArray<LatticeFermion>  > chrono_buf;
    const int N5;

  public:
    LinearExtrapolation5DChronoPredictor(const int N5_) : chrono_buf( new CircularBufferArray<LatticeFermion>(2, N5_) ), N5(N5_) {}

      
    ~LinearExtrapolation5DChronoPredictor(void) {}
    
    // Copying
    LinearExtrapolation5DChronoPredictor(const LinearExtrapolation5DChronoPredictor& p) : chrono_buf(p.chrono_buf), N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& A,
		    const multi1d<LatticeFermion>& chi) 
    { 
      START_CODE();

      switch( chrono_buf->size() ) { 
      case 0:
      {
	QDPIO::cout << "LinearExtrapolationPredictor: giving you zero" << endl;
	psi = zero;
      }
      break;

      case 1:
      { 
	QDPIO::cout << "LinearExtrapolationPredictor: giving you last soln" << endl;
	chrono_buf->get(0, psi);
      }
      break;

      case 2:
      {
	  
	QDPIO::cout << "LinearExtrapolationPredictor: giving you linear extrapolation" << endl;

	multi1d<LatticeFermion> y0(N5);
	chrono_buf->get(0, y0);

	multi1d<LatticeFermion> y1(N5);
	chrono_buf->get(1, y1);

	psi.resize(N5);
	for(int s = 0; s < N5; s++) { 
	  psi[s] = Real(2)*y0[s] - y1[s];         // Linear Extrapolation
	}
      }
      break;
      default:
	QDPIO::cerr << "Unknown case reached in LinearExtrapPredictor " << endl;
	QDP_abort(1);
	break;
      }

      END_CODE();
    }
    

    // No internal state so reset is a Nop
    void reset(void) {
      chrono_buf->reset();
    }

    // Ignore new vector
    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) 
    {
      START_CODE();

      QDPIO::cout << "LinearExtrapolationPredictor: registering new solution. ";
      chrono_buf->push(psi);
      QDPIO::cout << " number of vectors stored is = " << chrono_buf->size() << endl;

      END_CODE();
    }
    
  };
  
} // End Namespace Chroma

#endif 
