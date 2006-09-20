// -*- C++ -*-
// $Id: ovext_neuberger_strategy.h,v 3.1 2006-09-20 20:27:59 edwards Exp $
/*! \file
 *  \brief Ovext Neuberger rescale strategy
 */

#ifndef ovext_neuberger_strategy_h
#define ovext_neuberger_strategy_h

#include "chromabase.h"
#include "actions/ferm/fermacts/ovext_tuning_strategy.h"

namespace Chroma {
    
  /*! @ingroup fermacts */
  namespace OvExtNeubergerStrategyEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }
  
  
  //! Ovext Neuberger rescale strategy
  /*! @ingroup fermacts */
  class OvExtNeubergerStrategy : public AbsOvExtTuningStrategy 
  {
  public: 
    // Destructor is automatic
    ~OvExtNeubergerStrategy(void) {}

    OvExtNeubergerStrategy(void) {}
    
    // Strategy: Rightmost col and Lowest Row contain
    //
    //  sqrt( ( 1 - mu )/2 )
    //
    // which means beta should be set to 
    //
    //  (1-mu)/(2 resP[i] )
    
    void operator()( multi1d<Real>& beta,
		     const Real& coeffP,
		     const multi1d<Real>& resP,
		     const multi1d<Real>& resQ,
		     const Real& Mass) const {

      int size=resP.size();
      beta.resize(size);
      for(int i=0; i < size; i++) { 
	beta[i] = (Real(1)-Mass)/(Real(2)*resP[i]);
      }
    }
  private:
  };
  
}

#endif
