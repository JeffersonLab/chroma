// -*- C++ -*-
// $Id: ovext_const_div_by_resp_strategy.h,v 3.1 2006-09-20 20:27:58 edwards Exp $
/*! \file
 *  \brief Ovext rescale strategy
 */

#ifndef ovext_const_div_by_resp_h
#define ovext_const_div_by_resp_h

#include "chromabase.h"
#include "actions/ferm/fermacts/ovext_tuning_strategy.h"

namespace Chroma 
{
    
  /*! @ingroup fermacts */
  namespace OvExtConstDivByResPStrategyEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }
  
  
  //! Ovext rescale strategy
  /*! @ingroup fermacts */
  class OvExtConstDivByResPStrategy : public AbsOvExtTuningStrategy 
  {
  public: 
    // Destructor is automatic
    ~OvExtConstDivByResPStrategy(void) {}
    
    // Construct from program
    OvExtConstDivByResPStrategy(const Real& tuningConstant_ ): tuningConstant(tuningConstant_) {}
    
    
    void operator()( multi1d<Real>& beta,
		     const Real& coeffP,
		     const multi1d<Real>& resP,
		     const multi1d<Real>& resQ,
		     const Real& Mass) const {
      int size=resP.size();
      beta.resize(size);
      for(int i=0; i < size; i++) { 
	beta[i] = Real(1)/ (tuningConstant*resP[i]);
      }
    }
  private:
    const Real tuningConstant;
  };
  
};

#endif
