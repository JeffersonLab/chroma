// -*- C++ -*-
/*! \file
 * \brief Pick channel for QUDA Predictor
 *
 * Predictors for HMC
 */

#ifndef __quda_predictor_h__
#define __quda_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "quda.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"

namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace QUDA4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
    int getAndIncrGlobalQUDAChronoIndex();
  }

  //! Zero initial guess predictor
  /*! @ingroup predictor */
  class QUDA4DChronoPredictor : 
    public AbsTwoStepChronologicalPredictor4D<LatticeFermion> 
  {
  public:
    QUDA4DChronoPredictor(int max_chrono, QudaPrecisionType prec) :
    	_max_chrono(max_chrono), _prec(prec) {
      // Check Static Channel IDs:
      _X_index = QUDA4DChronoPredictorEnv::getAndIncrGlobalQUDAChronoIndex(); 
      _Y_index = QUDA4DChronoPredictorEnv::getAndIncrGlobalQUDAChronoIndex();
      QDPIO::cout << "Initing QUDAChrono Predictor with Channel IDs (X , Y)=("
		  << _X_index <<" , " << _Y_index << ") with max_chrono = "
		  << _max_chrono << std::endl;
    }

    // Destructor is automagic
    ~QUDA4DChronoPredictor(void) {}

    // Zero out psi -- it is a zero guess after all
    void predictX(LatticeFermion& X,
		  const LinearOperator<LatticeFermion> &A,
		  const LatticeFermion& chi) override
    {
	START_CODE();
	QDPIO::cout << "This is a special interface to use QUDA's predictor\n";
	QDPIO::cout << "predictX will return Zero\n";
	X = zero;
        END_CODE();
    }

    void predictY(LatticeFermion& Y,
		  const LinearOperator<LatticeFermion> &A,
		  const LatticeFermion& chi) override
    {
      START_CODE();
      QDPIO::cout << "This is a special interface to use QUDA's predictor\n";
      QDPIO::cout << "predictX will return Zero\n";
      Y = zero;
      END_CODE();
    }
    


    // No internal state so reset is a nop
    void reset(void) override {
      // Reset /flush an individual channel
      flushChronoQuda( getXIndex() );
      flushChronoQuda( getYIndex() );
    }

    // Ignore new std::vector
    void newXVector(const LatticeFermion& psi) override {
      QDPIO::cout << "QUDAPredictor: registering new X solution (not)" << std::endl;
    }

    void newYVector(const LatticeFermion& psi) override {
      QDPIO::cout << "QUDAPredictor: registering new Y solution (not)" << std::endl;
    }


    inline
    int getXIndex() const {
      return _X_index;
    }

    inline
    int  getYIndex() const {
      return _Y_index;
    }

    inline
    int getMaxChrono() const { 
      return _max_chrono;
    }
    
    inline
	QudaPrecisionType getChronoPrecision() const {
    	return _prec;
    }
  private:
    int _X_index;
    int _Y_index;
    int _max_chrono;
    QudaPrecisionType  _prec;
  };
  
} // End Namespace Chroma

#endif 
