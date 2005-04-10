// -*- C++ -*-
// $Id: linear_extrap_predictor.h,v 1.3 2005-04-10 21:46:43 edwards Exp $
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
  namespace LinearExtrapolation4DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  //! Last solution predictor
  /*! @ingroup predictor */
  class LinearExtrapolation4DChronoPredictor : 
    public AbsChronologicalPredictor4D<LatticeFermion> {

  private:
    Handle< CircularBuffer<LatticeFermion> > chrono_buf;

  public:

    LinearExtrapolation4DChronoPredictor(void) : chrono_buf(new CircularBuffer<LatticeFermion>((unsigned int)2)) {}

    // Destructor is automagic
    ~LinearExtrapolation4DChronoPredictor(void) {}


    void operator()(LatticeFermion& psi,
		    const LinearOperator<LatticeFermion>& A, 
		    const LatticeFermion& chi) {
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

    }
    
    // No internal state so reset is a nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
      chrono_buf->reset();
    }

    // Ignore new vector
    void newVector(const LatticeFermion& psi) {
      QDPIO::cout << "LinearExtrapolationPredictor: registering new solution. ";
      chrono_buf->push(psi);
      QDPIO::cout << " number of vectors stored is = " << chrono_buf->size() << endl;
    }

  };

  

  /*! @ingroup predictor */
  namespace LinearExtrapolation5DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };
  
  //! Last solution predictor
  /*! @ingroup predictor */
  class LinearExtrapolation5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {
    
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
		    const LinearOperator< multi1d<LatticeFermion> >& A,
		    const multi1d<LatticeFermion>& chi) { 

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


    }
    

    // No internal state so reset is a Nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
      chrono_buf->reset();

    }

    // Ignore new vector
    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) {
      QDPIO::cout << "LinearExtrapolationPredictor: registering new solution. ";
      chrono_buf->push(psi);
      QDPIO::cout << " number of vectors stored is = " << chrono_buf->size() << endl;
    }
    
  };
  
}; // End Namespace Chroma

#endif 
