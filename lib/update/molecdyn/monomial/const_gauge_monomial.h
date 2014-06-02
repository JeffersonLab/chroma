// -*- C++ -*-
/*! \file
 *  \brief Generic gauge action monomial wrapper
 */

#ifndef __const_gaugeact_monomial_h__
#define __const_gaugeact_monomial_h__

#include "chromabase.h"
#include "update/molecdyn/monomial/gauge_monomial.h"

namespace Chroma 
{
  /*! @ingroup monomial */
  namespace ConstGaugeMonomialEnv 
  {
    extern const string name;
    bool registerAll();
  }


  //! Wrapper class for  gauge monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class ConstGaugeMonomial : 
    public GaugeMonomial
  {
   public: 
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct out of a parameter struct. Check against the desired GaugeAct name
ConstGaugeMonomial(const GaugeMonomialParams& param_):GaugeMonomial(param_){}

    //! Create a suitable state and compute F
    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      GaugeMonomial::dsdq(F,s);
      ColorMatrix CF ;
      for(int mu(0);mu<Nd;mu++){
	CF=sum(F[mu])/toDouble(Layout::vol());F[mu]=CF;
      }
    }

    private:
      // Hide empty constructor and =
      ConstGaugeMonomial();
      void operator=(const ConstGaugeMonomial&);

    };


}; //end namespace chroma

#endif
