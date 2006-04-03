// $Id: ischiral_w.cc,v 3.0 2006-04-03 04:58:57 edwards Exp $
#include <chromabase.h>
#include "meas/eig/eig_w.h"

namespace Chroma {

Chirality isChiralVector(const LatticeFermion& chi)
{

  Chirality ret_val;

  Real chi_sq = Real(norm2(chi));

 
  //  LatticeFermion tmp = adj(chi)*Gamma(G5)*chi;
  Real chirality = real(sum(innerProduct(chi,Gamma(Ns*Ns-1)*chi,all)));

  // If chi is chiral, then 
  //
  //  chirality =  Sum Real Trace adj(chi)*Gamma_5*chi
  //   
  //            =  +/- Sum Real Trace adj(chi) chi
  //            =  +/- || chi ||^2   (+ for +ve chirality, - for -ve)
  //
  // So if  || chi ||^2 - | chirality | = 0 then 
  // we have a definite chirality.
  // 
  Real tmp1 = chi_sq - fabs(chirality);

  // To get a boolean out of < operator I have to apply
  // toBool. Is this because otherwise it is some kind of selector
  // for a mask?
  if ( toBool(tmp1 > fabs(chirality)*fuzz) ) {
    ret_val = CH_NONE;
  }
  else { 

    // To get a boolean out of < operator I have to apply
    // toBool. Is this because otherwise it is some kind of selector
    // for a mask?
    if( toBool( chirality > 0 ) ) { 
      ret_val = CH_PLUS;
    }
    else {
      ret_val = CH_MINUS;
    }
  }

  return ret_val;
}

}  // end namespace Chroma
