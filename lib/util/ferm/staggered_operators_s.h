// -*- C++ -*-
// $Id: staggered_operators_s.h,v 1.1 2006-11-18 05:41:51 kostas Exp $
/*! \file
 *  \brief Staggered  operators
 *
 */

#ifndef __staggered_operators_h__
#define __staggered_operators_h__


namespace Chroma 
{

  void SymShift(LatticeStaggeredPropagator& dest,const LatticeStaggeredPropagator& src,const  multi1d<LatticeColorMatrix>& u,const int mu) ;

  void SpinScalar(LatticeStaggeredPropagator& dest, const LatticeStaggeredPropagator& src, const  multi1d<LatticeColorMatrix>& u ) ;
  void SpinVector(LatticeStaggeredPropagator& dest, const LatticeStaggeredPropagator& src, const  multi1d<LatticeColorMatrix>& u,const int mu) ;
  void SpinTensor(LatticeStaggeredPropagator& dest, const LatticeStaggeredPropagator& src, const  multi1d<LatticeColorMatrix>& u,const int mu,const int nu) ;
  void SpinAxialVector (LatticeStaggeredPropagator& dest, const LatticeStaggeredPropagator& src, const  multi1d<LatticeColorMatrix>& u,const int mu) ;
  void SpinPseudoScalar(LatticeStaggeredPropagator& dest, const LatticeStaggeredPropagator& src, const  multi1d<LatticeColorMatrix>& u) ;

  void FlavorScalar(LatticeStaggeredPropagator& dest,const LatticeStaggeredPropagator& src,const  multi1d<LatticeColorMatrix>& u ) ;
  void FlavorVector(LatticeStaggeredPropagator& dest,const LatticeStaggeredPropagator& src,const  multi1d<LatticeColorMatrix>& u,const int mu) ;
  void FlavorTensor(LatticeStaggeredPropagator& dest,const LatticeStaggeredPropagator& src,const  multi1d<LatticeColorMatrix>& u,const int mu,const int nu) ;
  void FlavorAxialVector (LatticeStaggeredPropagator& dest,const LatticeStaggeredPropagator& src,const  multi1d<LatticeColorMatrix>& u,const int mu) ;
  void FlavorPseudoScalar(LatticeStaggeredPropagator& dest,const LatticeStaggeredPropagator& src,const  multi1d<LatticeColorMatrix>& u) ;

}  // end namespace Chroma

#endif
