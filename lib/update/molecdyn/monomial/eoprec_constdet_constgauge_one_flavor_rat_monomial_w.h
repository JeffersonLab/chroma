// -*- C++ -*-
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_constgauge_one_flavor_rat_monomial_w_h__
#define __prec_constgauge_one_flavor_rat_monomial_w_h__

#include "update/molecdyn/monomial/eoprec_constdet_one_flavor_rat_monomial_w.h"

namespace Chroma{

/*! @ingroup monomial */
namespace EvenOddPrecConstDetConstGaugeOneFlavorWilsonTypeFermRatMonomialEnv {
bool registerAll();
}


//! Wrapper class for  2-flavor even-odd prec ferm monomials
/*! @ingroup monomial
 *
 * Monomial is expected to be the same for these fermacts
 */
class EvenOddPrecConstDetConstGaugeOneFlavorWilsonTypeFermRatMonomial :
  public   EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial
{
   public: 
   // Typedefs to save typing
   typedef LatticeFermion               T;
   typedef multi1d<LatticeColorMatrix>  P;
   typedef multi1d<LatticeColorMatrix>  Q;

    // Construct out of a parameter struct. Check against the desired FermAct name
EvenOddPrecConstDetConstGaugeOneFlavorWilsonTypeFermRatMonomial(const OneFlavorWilsonTypeFermRatMonomialParams& param_): EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial(param_){}

virtual void dsdq(P& F, const AbsFieldState<P,Q>& s){
EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial::dsdq(F,s);

// Get at the ferion action 
 const WilsonTypeFermAct<T,P,Q>& FA = getFermAct();

 //FA.getFermBC().modify(F);
 ColorMatrix CF ;
 for(int mu(0);mu<Nd;mu++){
   CF=sum(F[mu])/Double( Layout::vol() ) ;
	 F[mu]=CF;
 }
} 

  private:
    // Hide empty constructor and =
    EvenOddPrecConstDetConstGaugeOneFlavorWilsonTypeFermRatMonomial();
void operator=(const EvenOddPrecConstDetConstGaugeOneFlavorWilsonTypeFermRatMonomial&);

  };


} //end namespace chroma

#endif
