// -*- C++ -*-
/*! \file
 *  \brief Ovext rescale strategy
 */

#ifndef ovext_constant_strategy_h
#define ovext_constant_strategy_h

#include "chromabase.h"
#include "actions/ferm/fermacts/ovext_tuning_strategy.h"

namespace Chroma 
{
    
  /*! @ingroup fermacts */
  namespace OvExtConstantStrategyEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }
  
  
  //! Ovext rescale strategy
  /*! @ingroup fermacts */
  class OvExtConstantStrategy : public AbsOvExtTuningStrategy 
  {
  public: 
    // Destructor is automatic
    ~OvExtConstantStrategy(void) {}
    
    // Construct from program
    OvExtConstantStrategy(const Real& tuningConstant_ ): tuningConstant(tuningConstant_) {}
    
    
    void operator()( multi1d<Real>& beta,
		     const Real& coeffP,
		     const multi1d<Real>& resP,
		     const multi1d<Real>& resQ,
		     const Real& Mass) const {
      QDPIO::cout << "Tuning constant is: " << tuningConstant << std::endl;
      int size=resP.size();
      beta.resize(size);
      for(int i=0; i < size; i++) { 
	beta[i] = tuningConstant;
	QDPIO::cout << "beta["<<i<<"] = " << beta[i] << std::endl;
      }
    }
  private:
    const Real tuningConstant;
  };
  
};

#endif
