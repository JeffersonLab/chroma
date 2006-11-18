// -*- C++ -*-
// $Id: staggered_operators_s.cc,v 1.1 2006-11-18 05:41:51 kostas Exp $
/*! \file
 *  \brief Staggered  operators
 *
 */

#include "chromabase.h"
#include "util/ferm/staggered_operators.h"

namespace Chroma 
{

  typedef LatticeStaggeredPropagator  T ;
  typedef multi1d<LatticeColorMatrix> G ;
  
  void SymShift(T& dest,const T& src,const G& u,const int mu){
    
  }

  void SpinScalar(T& dest, const T& src, const G& u ){

  }
  void SpinVector(T& dest, const T& src, const G& u,const int mu){

  }
  void SpinTensor(T& dest, const T& src, const G& u,const int mu,const int nu){

  }
  void SpinAxialVector (T& dest, const T& src, const G& u,const int mu){

  }
  void SpinPseudoScalar(T& dest, const T& src, const G& u){

  }

  void FlavorScalar(T& dest,const T& src,const G& u ){

  }
  
  void FlavorVector(T& dest,const T& src,const G& u,const int mu){

  }
  
  void FlavorTensor(T& dest,const T& src,const G& u,const int mu,const int nu){

  }
  
  void FlavorAxialVector (T& dest,const T& src,const G& u,const int mu){

  }

  void FlavorPseudoScalar(T& dest,const T& src,const G& u){

  }

}  // end namespace Chroma

#endif
