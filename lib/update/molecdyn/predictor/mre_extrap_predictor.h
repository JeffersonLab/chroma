#ifndef __mre_extrap_predictor_h__
#define __mre_extrap_predictor_h__

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/circular_buffer.h"

namespace Chroma { 
  
  namespace MinimalResidualExtrapolation4DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  class MinimalResidualExtrapolation4DChronoPredictor  
    : public AbsChronologicalPredictor4D<LatticeFermion> {

  private:
    Handle< CircularBuffer<LatticeFermion> > chrono_buf;
    
    void find_extrap_solution(LatticeFermion& psi, 
			      const LinearOperator<LatticeFermion>& A,
			      const LatticeFermion& chi);
  public:
    
    MinimalResidualExtrapolation4DChronoPredictor(unsigned int max_chrono) : chrono_buf(new CircularBuffer<LatticeFermion>(max_chrono)) {}
    
    // Destructor is automagic
    ~MinimalResidualExtrapolation4DChronoPredictor(void) {}

    // Do the hard work
    void operator()(LatticeFermion& psi, 
		    const LinearOperator<LatticeFermion>& A,
		    const LatticeFermion& chi);
    
    // No internal state so reset is a nop
    void reset(void) {
      QDPIO::cout << "Resetting Chrono Predictor" << endl;
      chrono_buf->reset();
    }

    // Ignore new vector
    void newVector(const LatticeFermion& psi) {
      QDPIO::cout << "MinimalResidualExtrapolationPredictor: registering new solution. ";
      chrono_buf->push(psi);
      QDPIO::cout << " number of vectors stored is = " << chrono_buf->size() << endl;
    }

  };

  
#if 0
  namespace MinimalResidualExtrapolation5DChronoPredictorEnv {
    extern const std::string name;
    extern const bool registered;
  };
  
  class MinimalResidualExtrapolation5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {
    
    private: 
    Handle< CircularBufferArray<LatticeFermion>  > chrono_buf;
    const int N5;

    public:
    MinimalResidualExtrapolation5DChronoPredictor(const int N5_) : chrono_buf( new CircularBufferArray<LatticeFermion>(2, N5_) ), N5(N5_) {}

      
    ~MinimalResidualExtrapolation5DChronoPredictor(void) {}
    
    // Copying
    MinimalResidualExtrapolation5DChronoPredictor(const MinimalResidualExtrapolation5DChronoPredictor& p) : chrono_buf(p.chrono_buf), N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi) { 

      switch( chrono_buf->size() ) { 
      case 0:
	{
	  QDPIO::cout << "MinimalResidualExtrapolationPredictor: giving you zero" << endl;
	  psi = zero;
	}
	break;

      case 1:
	{ 
	  QDPIO::cout << "MinimalResidualExtrapolationPredictor: giving you last soln" << endl;
	  chrono_buf->get(0, psi);
	}
	break;

      case 2:
	{
	  
	  QDPIO::cout << "MinimalResidualExtrapolationPredictor: giving you linear extrapolation" << endl;

	  multi1d<LatticeFermion> y0(N5);
	  chrono_buf->get(0, y0);

	  multi1d<LatticeFermion> y1(N5);
	  chrono_buf->get(1, y1);

	  psi.resize(N5);
	  for(int s = 0; s < N5; s++) { 
	    psi[s] = Real(2)*y0[s] - y1[s];         // MinimalResidual Extrapolation
	  }
	}
	break;
      default:
	QDPIO::cerr << "Unknown case reached in MinimalResidualExtrapPredictor " << endl;
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
      QDPIO::cout << "MinimalResidualExtrapolationPredictor: registering new solution. ";
      chrono_buf->push(psi);
      QDPIO::cout << " number of vectors stored is = " << chrono_buf->size() << endl;
    }


  };
#endif
  
}; // End Namespace Chroma

#endif 
